const _libcapgen = window.libcapgen();

let instance = null;
const callbacks = [];
_libcapgen.onRuntimeInitialized = () => {
  instance = {
    optimize: (...args) => _libcapgen.optimize(...args)
  };
  callbacks.map(c => c());
};

export const libcapgen = () =>
  new Promise(resolve => {
    if (instance === null) {
      callbacks.push(() => resolve(instance));
    } else {
      resolve(instance);
    }
  });
