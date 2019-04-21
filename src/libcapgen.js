const _libcapgen = window.libcapgen();

let instance = null;
const callbacks = [];
_libcapgen.onRuntimeInitialized = () => {
  instance = _libcapgen;
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
